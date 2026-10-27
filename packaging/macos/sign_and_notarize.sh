#!/bin/bash
# Sign the built app with your Developer ID, have Apple notarize it, staple the ticket, then make and notarize the
# DMG, so that the app opens on any Mac without "unidentified developer" warnings.
# Once per Mac:
#   1. a "Developer ID Application" certificate in the login keychain (Xcode > Settings > Accounts > Manage
#      Certificates > + > Developer ID Application; only the account holder can create one). Check with
#        security find-identity -v -p codesigning
#   2. the notary credentials in the keychain (an app-specific password from appleid.apple.com > Sign-In and
#      Security > App-Specific Passwords; the Team ID is under developer.apple.com/account > Membership):
#        xcrun notarytool store-credentials stannotator --apple-id YOU@EXAMPLE.COM --team-id TEAMID
# Then, after packaging/macos/build_app.sh (or run build_app.sh with DEVELOPER_ID set, which calls this):
#   DEVELOPER_ID="Developer ID Application: Your Name (TEAMID)" bash packaging/macos/sign_and_notarize.sh
set -euo pipefail
cd "$(dirname "$0")/../.."
APP="dist/STAnnotator.app"
ENT=packaging/macos/entitlements.plist
PROFILE=${NOTARY_PROFILE:-stannotator}
: "${DEVELOPER_ID:?set DEVELOPER_ID to your 'Developer ID Application: ...' identity (security find-identity -v -p codesigning)}"
[ -d "$APP" ] || { echo "no $APP: run packaging/macos/build_app.sh first"; exit 1; }
case "$DEVELOPER_ID" in
  "Developer ID Application:"*) ;;
  *) echo "DEVELOPER_ID is \"$DEVELOPER_ID\": Apple notarizes only apps signed with a \"Developer ID Application\""
     echo "certificate (\"Apple Development\" and \"Apple Distribution\" are rejected). Create one in Xcode > Settings >"
     echo "Accounts > Manage Certificates > + > Developer ID Application, then use its whole name."; exit 1 ;;
esac
if ! security find-identity -v -p codesigning | grep -qF "$DEVELOPER_ID"; then
  echo "the identity \"$DEVELOPER_ID\" is not in your keychain; these are:"; security find-identity -v -p codesigning
  if ! security find-identity -v -p codesigning | grep -q "Developer ID Application:"; then
    echo "None is a \"Developer ID Application\" certificate, the kind that signs apps for other Macs (\"Apple"
    echo "Development\" and \"Apple Distribution\" are for Xcode and the App Store). Create it in Xcode > Settings >"
    echo "Accounts > Manage Certificates > + > Developer ID Application (account holder only)."
  else
    echo "Set DEVELOPER_ID to the whole quoted name of the \"Developer ID Application: ...\" line above."
  fi
  exit 1
fi
if ! xcrun notarytool history --keychain-profile "$PROFILE" >/dev/null 2>&1; then
  echo "no notary credentials named $PROFILE: run once"
  echo "  xcrun notarytool store-credentials $PROFILE --apple-id YOU@EXAMPLE.COM --team-id TEAMID"; exit 1
fi
notary_failed() {  # Apple's reasons, file by file
  id=$(grep -m1 -o "id: [0-9a-f-]*" dist/notarize.log | cut -d" " -f2)
  echo "Apple rejected it. Its log ($id):"
  xcrun notarytool log "$id" --keychain-profile "$PROFILE" | tee dist/notary_rejection.json \
    | grep -E '"(path|message)"' | head -40
  echo "(the whole log: dist/notary_rejection.json)"
  exit 1
}
echo "signing $APP"
# sign every bundled library and executable from the inside out, then the app itself
find "$APP/Contents" -type f -print0 | while IFS= read -r -d '' f; do
  if file -b "$f" | grep -q "Mach-O"; then
    codesign --force --timestamp --options runtime --entitlements "$ENT" --sign "$DEVELOPER_ID" "$f"
  fi
done
codesign --force --timestamp --options runtime --entitlements "$ENT" --sign "$DEVELOPER_ID" "$APP"
codesign --verify --deep --strict --verbose=2 "$APP"
ditto -c -k --keepParent "$APP" dist/notarize.zip
echo "notarizing the app (Apple usually answers within minutes)"
xcrun notarytool submit dist/notarize.zip --keychain-profile "$PROFILE" --wait | tee dist/notarize.log
grep -q "status: Accepted" dist/notarize.log || notary_failed
xcrun stapler staple "$APP"
rm -f dist/notarize.zip dist/STAnnotator.dmg
hdiutil create -volname "STAnnotator" -srcfolder "$APP" -ov -format ULMO dist/STAnnotator.dmg \
  || hdiutil create -volname "STAnnotator" -srcfolder "$APP" -ov -format UDZO dist/STAnnotator.dmg
codesign --force --timestamp --sign "$DEVELOPER_ID" dist/STAnnotator.dmg
echo "notarizing the DMG"
xcrun notarytool submit dist/STAnnotator.dmg --keychain-profile "$PROFILE" --wait | tee dist/notarize.log
grep -q "status: Accepted" dist/notarize.log || notary_failed
rm -f dist/notarize.log
xcrun stapler staple dist/STAnnotator.dmg
spctl --assess --type open --context context:primary-signature -v dist/STAnnotator.dmg
echo "signed and notarized: dist/STAnnotator.dmg"
