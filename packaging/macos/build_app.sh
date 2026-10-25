#!/bin/bash
# Build the standalone app on a Mac with Apple silicon: dist/Spatial Annotator.app and dist/SpatialAnnotator.dmg.
# Signing and notarization (optional, needs an Apple Developer ID): packaging/macos/sign_and_notarize.sh
set -euo pipefail
cd "$(dirname "$0")/../.."
[ "$(uname -m)" = "arm64" ] || { echo "build on an Apple silicon Mac"; exit 1; }
uv venv --allow-existing .build-venv --python 3.12
# --reinstall-package: an earlier build left the same version in .build-venv; always take this folder's code
uv pip install --python .build-venv/bin/python --reinstall-package celltypeagent ".[app]" pyinstaller
# app icon (.icns) from the PNG
ICONSET=$(mktemp -d)/AppIcon.iconset; mkdir -p "$ICONSET"
for s in 16 32 128 256 512; do
  sips -z $s $s celltypeagent/data/app_icon.png --out "$ICONSET/icon_${s}x${s}.png" >/dev/null
  sips -z $((2*s)) $((2*s)) celltypeagent/data/app_icon.png --out "$ICONSET/icon_${s}x${s}@2x.png" >/dev/null
done
iconutil -c icns "$ICONSET" -o packaging/macos/AppIcon.icns
.build-venv/bin/pyinstaller --noconfirm --clean packaging/spatial_annotator.spec
codesign --force --deep --sign - "dist/Spatial Annotator.app"   # ad-hoc signature (replaced when notarizing)
rm -f dist/SpatialAnnotator.dmg
hdiutil create -volname "Spatial Annotator" -srcfolder "dist/Spatial Annotator.app" -ov -format UDZO dist/SpatialAnnotator.dmg
du -sh "dist/Spatial Annotator.app" dist/SpatialAnnotator.dmg
echo "built: $(cat celltypeagent/data/build.txt)"
# a launcher made by 'celltypeagent make-app' has the same name and runs an older source folder
if [ -d "$HOME/Applications/Spatial Annotator.app" ]; then
  echo "Note: $HOME/Applications/Spatial Annotator.app is a launcher of a source install ($(grep -o '/[^\"]*/.venv' "$HOME/Applications/Spatial Annotator.app/Contents/MacOS/"* 2>/dev/null | head -1))."
  echo "      Spotlight and the Dock may open it instead of this build: move it to the Trash to avoid mixing them up."
fi
