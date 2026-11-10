# PyInstaller spec of the standalone STAnnotator app.
#   macOS (Apple silicon): packaging/macos/build_app.sh  ->  dist/STAnnotator.app
#   Windows 10 / 11 (x64): packaging\windows\build.ps1  ->  dist\STAnnotator\STAnnotator.exe (the window) and
#                          stagent.exe (console: Claude Code runs "stagent.exe mcp"), zipped as dist\STAnnotator-windows.zip
#   Linux (to test the bundle): pyinstaller packaging/stannotator.spec  ->  dist/STAnnotator/
import sys
from pathlib import Path
from PyInstaller.utils.hooks import collect_all, copy_metadata

ROOT = Path(SPECPATH).parent
VERSION = __import__("re").search(r'__version__ = "([^"]+)"', (ROOT / "stagent" / "__init__.py").read_text()).group(1)
if not (ROOT / "stagent" / "data" / "report_upload_url.txt").exists():
    print("NOTE: no stagent/data/report_upload_url.txt: Report bug will ask testers for the upload link")
datas = [(str(ROOT / "stagent" / "static"), "stagent/static"),
         (str(ROOT / "stagent" / "data"), "stagent/data")]
# source of the modules an exported notebook carries (the bundle otherwise holds only compiled code)
import ast, re
# the list may run over several lines; UTF-8, not the Windows code page
CARRIED = ast.literal_eval(re.search(r"^CARRIED = (\[.*?\])", (ROOT / "stagent" / "notebook.py").read_text(encoding="utf-8"),
                                     re.M | re.S).group(1))
datas += [(str(ROOT / "stagent" / f"{m}.py"), "stagent/algorithm_src") for m in CARRIED]
# backend_pdf: Save figure writes PDF; PyInstaller's matplotlib hook packs only the backends named in matplotlib.use()
binaries, hidden = [], ["stagent.mcp_server", "stagent.webapp", "stagent.desktop", "matplotlib.backends.backend_pdf"]

# Left out to keep the app small; none of it is ever imported by the app, so nothing it does changes:
# - imagecodecs: the codecs of other file formats (CZI, FITS, AVIF, ...). Images are read with tifffile, which can call
#   only the TIFF compressions below (jpeg8 and ljpeg through imagecodecs.jpeg_decode); imagecodecs loads a codec
#   module when it is first used.
TIFF_CODECS = {"_ccitt", "_deflate", "_imcd", "_jetraw", "_jpeg2k", "_jpeg8", "_jpegxl", "_jpegxr", "_lerc", "_ljpeg",
               "_lzma", "_png", "_webp", "_zlib", "_zstd", "_shared", "_shared_cython"}
# - pyarrow: Flight (RPC) and the S3 file system (optional in pyarrow.fs); the app reads local parquet files only.
# - cryptography (pulled in by optional parts of mcp), numba's CUDA target (no NVIDIA GPU on a Mac), test suites,
#   C / Cython sources and headers.
DROP_MODULES = ["cryptography", "pyarrow.flight", "pyarrow._flight", "pyarrow._s3fs", "numba.cuda"]
DROP_FILES = re.compile(r"(^|[/\\])(tests|include)([/\\]|$)|\.(pyx|pxd|pxi|h|hpp|c|cpp)$|libarrow_flight|libarrow_s3|"
                        r"pyarrow[/\\]_(flight|s3fs)\.|^cryptography|[/\\]cryptography[/\\]|numba[/\\]cuda[/\\]")


def dropped_module(name):
    if name.startswith("imagecodecs._") and name.split(".")[1] not in TIFF_CODECS:
        return True
    return ".tests" in name or any(name == m or name.startswith(m + ".") for m in DROP_MODULES)


def dropped_file(dest):
    p = Path(dest)
    if p.parent.name == "imagecodecs" and p.name.startswith("_") and p.name.split(".")[0] not in TIFF_CODECS | {"__init__", "__main__"} \
            and p.suffix in (".so", ".pyd", ".dylib"):
        return True
    return bool(DROP_FILES.search(str(dest)))


# packages that load parts of themselves dynamically (codecs, event loops, numba-compiled modules, plugins)
# (Windows: pywebview drives Edge WebView2 through pythonnet / clr_loader, which load .NET assemblies and DLLs)
for pkg in ["zarr", "numcodecs", "imagecodecs", "tifffile", "pynndescent", "numba", "llvmlite", "igraph", "umap", "sklearn", "harmonypy",
            "uvicorn", "mcp", "webview", "pydeseq2", "anndata", "formulaic", "formulaic_contrasts"] + (["clr_loader", "pythonnet"] if sys.platform == "win32" else []):
    try:
        d, b, h = collect_all(pkg)
    except Exception:
        continue
    datas += [x for x in d if not dropped_file(x[1] + "/" + Path(x[0]).name)]
    binaries += [x for x in b if not dropped_file(x[1] + "/" + Path(x[0]).name)]
    hidden += [m for m in h if not dropped_module(m)]
for pkg in ["zarr", "numcodecs", "mcp", "imagecodecs", "tifffile", "uvicorn", "starlette", "pydantic"]:
    try:
        datas += copy_metadata(pkg)
    except Exception:
        pass

a = Analysis([str(ROOT / "packaging" / "entry.py")], pathex=[str(ROOT)], datas=datas, binaries=binaries,
             hiddenimports=hidden, excludes=["torch", "tkinter", "IPython", "pytest", "notebook", "tensorflow", "keras"] + DROP_MODULES
             + [f"imagecodecs.{p.stem.split('.')[0]}" for p in Path(__import__("imagecodecs").__file__).parent.glob("_*")
                if p.suffix in (".so", ".pyd") and p.stem.split(".")[0] not in TIFF_CODECS], noarchive=False)
a.datas = [x for x in a.datas if not dropped_file(x[0])]
a.binaries = [x for x in a.binaries if not dropped_file(x[0])]
pyz = PYZ(a.pure)
if sys.platform == "win32":
    # two programs over one bundle: the window has no console; Claude Code needs a console program for the MCP server
    # (stdin / stdout), which config.self_command finds beside the window's
    icon = str(ROOT / "stagent" / "data" / "app_icon.png")  # converted to .ico by PyInstaller (Pillow)
    exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="STAnnotator", console=False, icon=icon)
    cli = EXE(pyz, a.scripts, [], exclude_binaries=True, name="stagent", console=True, icon=icon)
    coll = COLLECT(exe, cli, a.binaries, a.datas, name="STAnnotator")
else:
    exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="STAnnotator", console=sys.platform != "darwin",
              argv_emulation=False, target_arch="arm64" if sys.platform == "darwin" else None)
    coll = COLLECT(exe, a.binaries, a.datas, name="STAnnotator")
if sys.platform == "darwin":
    icon = ROOT / "packaging" / "macos" / "AppIcon.icns"
    app = BUNDLE(coll, name="STAnnotator.app", icon=str(icon) if icon.exists() else None,
                 bundle_identifier="io.github.caokai1073.stannotator",
                 info_plist={"CFBundleName": "STAnnotator", "CFBundleDisplayName": "STAnnotator",
                             "CFBundleShortVersionString": VERSION, "LSMinimumSystemVersion": "11.0",
                             "NSHighResolutionCapable": True})
