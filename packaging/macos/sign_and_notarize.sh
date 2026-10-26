#!/bin/bash
# Sign the built app with a Developer ID, notarize it with Apple and staple the ticket, then rebuild the DMG.
# Needs an Apple Developer account. Once per Mac, store the notary credentials in the keychain:
#   xcrun notarytool store-credentials spatial-annotator --apple-id YOU@EXAMPLE.COM --team-id TEAMID
#     (it asks for an app-specific password from appleid.apple.com)
# Then, after packaging/macos/build_app.sh:
#        DEVELOPER_ID="Developer ID Application: Your Name (TEAMID)" packaging/macos/sign_and_notarize.sh
set -euo pipefail
cd "$(dirname "$0")/../.."
APP="dist/Spatial Annotator.app"
ENT=packaging/macos/entitlements.plist
PROFILE=${NOTARY_PROFILE:-spatial-annotator}
: "${DEVELOPER_ID:?set DEVELOPER_ID to your 'Developer ID Application: ...' identity}"
# sign every bundled library and executable from the inside out, then the app itself
find "$APP/Contents" -type f -print0 | while IFS= read -r -d '' f; do
  if file -b "$f" | grep -q "Mach-O"; then
    codesign --force --timestamp --options runtime --entitlements "$ENT" --sign "$DEVELOPER_ID" "$f"
  fi
done
codesign --force --timestamp --options runtime --entitlements "$ENT" --sign "$DEVELOPER_ID" "$APP"
codesign --verify --deep --strict --verbose=2 "$APP"
ditto -c -k --keepParent "$APP" dist/notarize.zip
xcrun notarytool submit dist/notarize.zip --keychain-profile "$PROFILE" --wait
xcrun stapler staple "$APP"
rm -f dist/notarize.zip dist/SpatialAnnotator.dmg
hdiutil create -volname "Spatial Annotator" -srcfolder "$APP" -ov -format ULMO dist/SpatialAnnotator.dmg \
  || hdiutil create -volname "Spatial Annotator" -srcfolder "$APP" -ov -format UDZO dist/SpatialAnnotator.dmg
codesign --force --timestamp --sign "$DEVELOPER_ID" dist/SpatialAnnotator.dmg
xcrun notarytool submit dist/SpatialAnnotator.dmg --keychain-profile "$PROFILE" --wait
xcrun stapler staple dist/SpatialAnnotator.dmg
spctl --assess --type open --context context:primary-signature -v dist/SpatialAnnotator.dmg
echo "signed and notarized: dist/SpatialAnnotator.dmg"
