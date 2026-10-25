# PyInstaller spec of the standalone Spatial Annotator app.
#   macOS (Apple silicon): packaging/macos/build_app.sh  ->  dist/Spatial Annotator.app
#   Linux (to test the bundle): pyinstaller packaging/spatial_annotator.spec  ->  dist/Spatial Annotator/
import sys
from pathlib import Path
from PyInstaller.utils.hooks import collect_all, copy_metadata

ROOT = Path(SPECPATH).parent
datas = [(str(ROOT / "celltypeagent" / "static"), "celltypeagent/static"),
         (str(ROOT / "celltypeagent" / "data"), "celltypeagent/data")]
binaries, hidden = [], ["celltypeagent.mcp_server", "celltypeagent.webapp", "celltypeagent.desktop"]
# packages that load parts of themselves dynamically (codecs, event loops, numba-compiled modules, plugins)
for pkg in ["zarr", "numcodecs", "imagecodecs", "tifffile", "pynndescent", "numba", "llvmlite", "igraph", "umap", "sklearn",
            "uvicorn", "mcp", "webview"]:
    try:
        d, b, h = collect_all(pkg)
    except Exception:
        continue
    datas += d; binaries += b; hidden += h
for pkg in ["zarr", "numcodecs", "mcp", "imagecodecs", "tifffile", "uvicorn", "starlette", "pydantic"]:
    try:
        datas += copy_metadata(pkg)
    except Exception:
        pass

a = Analysis([str(ROOT / "packaging" / "entry.py")], pathex=[str(ROOT)], datas=datas, binaries=binaries,
             hiddenimports=hidden, excludes=["torch", "tkinter", "IPython", "pytest", "notebook", "tensorflow", "keras"], noarchive=False)
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="Spatial Annotator", console=sys.platform != "darwin",
          argv_emulation=False, target_arch="arm64" if sys.platform == "darwin" else None)
coll = COLLECT(exe, a.binaries, a.datas, name="Spatial Annotator")
if sys.platform == "darwin":
    icon = ROOT / "packaging" / "macos" / "AppIcon.icns"
    app = BUNDLE(coll, name="Spatial Annotator.app", icon=str(icon) if icon.exists() else None,
                 bundle_identifier="io.github.caokai1073.spatialannotator",
                 info_plist={"CFBundleName": "Spatial Annotator", "CFBundleDisplayName": "Spatial Annotator",
                             "CFBundleShortVersionString": "0.2.0", "LSMinimumSystemVersion": "11.0",
                             "NSHighResolutionCapable": True})
