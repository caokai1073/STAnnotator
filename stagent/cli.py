"""Command line: stagent {app, make-app, import, import-xenium, import-visium, samples, viewer, mcp}."""
import argparse, os, platform, shutil, sys
from pathlib import Path
from . import config


APP_PLIST = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0"><dict>
  <key>CFBundleName</key><string>STAnnotator</string>
  <key>CFBundleDisplayName</key><string>STAnnotator</string>
  <key>CFBundleIdentifier</key><string>io.github.caokai1073.stannotator</string>
  <key>CFBundleExecutable</key><string>STAnnotator</string>
  <key>CFBundlePackageType</key><string>APPL</string>
  <key>CFBundleShortVersionString</key><string>{version}</string>
  <key>CFBundleIconFile</key><string>AppIcon</string>
  <key>LSMinimumSystemVersion</key><string>11.0</string>
  <key>NSHighResolutionCapable</key><true/>
</dict></plist>
"""


def make_app(dest=None):
    """Build a small macOS app bundle that starts this installation's desktop app (no copy of Python inside: it
    runs the environment it was made from, so keep that folder where it is)."""
    import subprocess, tempfile
    if platform.system() != "Darwin":
        raise SystemExit("make-app builds a macOS app; on this system run: stagent app")
    dest = Path(dest or Path.home() / "Applications").expanduser()
    app = dest / "STAnnotator.app"
    if app.exists():
        shutil.rmtree(app)
    (app / "Contents" / "MacOS").mkdir(parents=True)
    (app / "Contents" / "Resources").mkdir()
    (app / "Contents" / "Info.plist").write_text(APP_PLIST.format(version=__import__("stagent").__version__))
    exe = app / "Contents" / "MacOS" / "STAnnotator"
    exe.write_text("#!/bin/sh\n"
                   f'export STAGENT_HOME="{config.HOME}"\n'
                   f'exec "{sys.executable}" -m stagent.desktop "$@" >> "{config.HOME / "app.log"}" 2>&1\n')
    exe.chmod(0o755)
    icon = config.PACKAGE / "data" / "app_icon.png"
    try:  # .icns from the PNG with the system's iconutil
        from PIL import Image
        with tempfile.TemporaryDirectory() as tmp:
            iconset = Path(tmp) / "AppIcon.iconset"
            iconset.mkdir()
            im = Image.open(icon)
            for s in (16, 32, 128, 256, 512):
                im.resize((s, s), Image.LANCZOS).save(iconset / f"icon_{s}x{s}.png")
                im.resize((2 * s, 2 * s), Image.LANCZOS).save(iconset / f"icon_{s}x{s}@2x.png")
            subprocess.run(["iconutil", "-c", "icns", str(iconset), "-o", str(app / "Contents" / "Resources" / "AppIcon.icns")],
                           check=True)
    except Exception as e:
        print(f"(no app icon: {e})")
    subprocess.run(["touch", str(app)], check=False)
    return app


def main(argv=None):
    ap = argparse.ArgumentParser(prog="stagent", description="Spatial cell type annotation with Claude.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    im = sub.add_parser("import", help="import one sample from h5ad / RDS + registered H&E + registration file")
    im.add_argument("--name", required=True, help="sample id used in the viewer and in Claude")
    im.add_argument("--h5ad", required=True, nargs="+", help="one or more h5ad files, or R files (.rds / .rda) with a Seurat or SpatialExperiment object (merged)")
    im.add_argument("--he", help="H&E image (pyramidal TIFF / OME-TIFF); without --registration it waits for an alignment")
    im.add_argument("--registration", help="affine CSV (a11,a12,a21,a22,tx,ty) or Xenium Explorer 3x3 matrix")
    im.add_argument("--section-column", help="obs column that identifies the section (if the h5ad holds several)")
    im.add_argument("--section", help="value of --section-column to import")
    im.add_argument("--x-column"); im.add_argument("--y-column")
    im.add_argument("--counts-layer", help="layer with raw counts (default: 'counts' if present, else X)")
    im.add_argument("--region-column", help="optional obs column with region annotations to start from")
    im.add_argument("--xenium-px-um", type=float, default=0.2125, help="morphology pixel size for Xenium Explorer matrices")
    im.add_argument("--description", default="")
    im.add_argument("--overwrite", action="store_true")
    ix = sub.add_parser("import-xenium", help="import one sample from a Xenium output folder, or only its files")
    ix.add_argument("--name", required=True, help="sample id used in the viewer and in Claude")
    ix.add_argument("--xenium", required=True, nargs="+", help="Xenium output folder, or its files: cells.parquet (or "
                    "cells.csv.gz) and cell_feature_matrix.h5, optionally experiment.xenium and the DAPI ome.tif")
    ix.add_argument("--he", help="optional H&E image registered with Xenium Explorer")
    ix.add_argument("--registration", help="the H&E's alignment file (Xenium Explorer) or an affine CSV")
    ix.add_argument("--description", default="")
    ix.add_argument("--overwrite", action="store_true")
    iv = sub.add_parser("import-visium", help="import a Visium HD (or Visium) Space Ranger output: bins analysed at 8 um")
    iv.add_argument("--name", required=True, help="sample id used in the viewer and in Claude")
    iv.add_argument("--path", required=True, help="outs/, binned_outputs/ or one square_XXXum folder (2 um bins are "
                    "pooled into 8 um; 8 um or more are used as they are; the whole output: its 8 um bins)")
    iv.add_argument("--he", help="optional full-resolution image Space Ranger was given (else tissue_hires_image.png)")
    iv.add_argument("--description", default="")
    iv.add_argument("--overwrite", action="store_true")
    sub.add_parser("samples", help="list imported samples")
    ap_ = sub.add_parser("app", help="start the desktop app (native window)")
    ap_.add_argument("--browser", action="store_true", help="use the default browser instead of a window")
    ma = sub.add_parser("make-app", help="macOS: create 'STAnnotator.app' in ~/Applications")
    ma.add_argument("--dest", help="folder for the app (default ~/Applications)")
    vw = sub.add_parser("viewer", help="start the web viewer")
    vw.add_argument("--port", type=int, default=8765)
    sub.add_parser("mcp", help="run the MCP server on stdio (Claude Code starts this itself)")
    rs = sub.add_parser("remote-serve", help="on a server: start the viewer in the background for the app on another "
                                             "computer (it connects over SSH) and print its port")
    rs.add_argument("--check", action="store_true", help="only print the version and data folder")
    a = ap.parse_args(argv)
    config.ensure_dirs()
    if a.cmd == "import":
        from .importer import import_sample
        import_sample(a.name, a.h5ad, a.he, a.registration, a.section_column, a.section, a.x_column, a.y_column,
                      a.counts_layer, a.region_column, a.xenium_px_um, a.description, a.overwrite)
    elif a.cmd == "import-xenium":
        from .importer import import_xenium
        import_xenium(a.name, a.xenium, a.he, a.registration, a.description, a.overwrite)
    elif a.cmd == "import-visium":
        from .importer import import_visium
        import_visium(a.name, a.path, a.he, a.description, a.overwrite)
    elif a.cmd == "app":
        from .desktop import main as app
        app(["--browser"] if a.browser else [])
    elif a.cmd == "make-app":
        print(f"created {make_app(a.dest)}\nOpen it from Finder or drag it to the Dock. Data folder: {config.HOME}")
    elif a.cmd == "samples":
        from .session import list_samples
        t = list_samples()
        print(t.to_string(index=False) if len(t) else f"no samples in {config.SAMPLES}")
    elif a.cmd == "viewer":
        sys.argv = ["stagent-web", "--port", str(a.port)]
        from .webapp import main as web
        web()
    elif a.cmd == "remote-serve":
        from .remote import serve
        serve(a.check)
    elif a.cmd == "mcp":
        from .mcp_server import main as mcp
        mcp()
        sys.stdout = open(os.devnull, "w")  # the client has closed stdout: nothing may flush it at exit


if __name__ == "__main__":
    main()
