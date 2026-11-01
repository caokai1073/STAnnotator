#!/bin/bash
# Sign the built app with your Developer ID, put it in a signed DMG, have Apple notarize the DMG (and so the app in
# it) and staple the ticket, so that the app opens on any Mac without "unidentified developer" warnings. A
# notarization still running when you stop the script is picked up with NOTARY_ID=<id> (no new upload).
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
[ -n "${NOTARY_ID:-}" ] || : "${DEVELOPER_ID:?set DEVELOPER_ID to your 'Developer ID Application: ...' identity (security find-identity -v -p codesigning)}"
DEVELOPER_ID=${DEVELOPER_ID:-}
[ -d "$APP" ] || { echo "no $APP: run packaging/macos/build_app.sh first"; exit 1; }
case "${NOTARY_ID:+resume}$DEVELOPER_ID" in
  resume*|"Developer ID Application:"*) ;;
  *) echo "DEVELOPER_ID is \"$DEVELOPER_ID\": Apple notarizes only apps signed with a \"Developer ID Application\""
     echo "certificate (\"Apple Development\" and \"Apple Distribution\" are rejected). Create one in Xcode > Settings >"
     echo "Accounts > Manage Certificates > + > Developer ID Application, then use its whole name."; exit 1 ;;
esac
if [ -z "${NOTARY_ID:-}" ] && ! security find-identity -v -p codesigning | grep -qF "$DEVELOPER_ID"; then
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
  echo "Apple rejected it. Its log ($1):"
  xcrun notarytool log "$1" --keychain-profile "$PROFILE" | tee dist/notary_rejection.json \
    | grep -E '"(path|message)"' | head -40
  echo "(the whole log: dist/notary_rejection.json)"
  exit 1
}
sign() {  # Apple's timestamp server sometimes does not answer one of the hundreds of requests ("A timestamp was
          # expected but was not found"): try each file up to 5 times
  local i
  for i in 1 2 3 4 5; do
    codesign --force --timestamp "$@" && return 0
    echo "retrying in $((i * 5)) s (attempt $i failed: usually Apple's timestamp server not answering)"
    sleep $((i * 5))
  done
  return 1
}
DMG=dist/STAnnotator.dmg
if [ -z "${NOTARY_ID:-}" ]; then
  echo "signing $APP"
  # sign every bundled library and executable from the inside out, then the app itself
  find "$APP/Contents" -type f -print0 | while IFS= read -r -d '' f; do
    if file -b "$f" | grep -q "Mach-O"; then
      sign --options runtime --entitlements "$ENT" --sign "$DEVELOPER_ID" "$f"
    fi
  done
  sign --options runtime --entitlements "$ENT" --sign "$DEVELOPER_ID" "$APP"
  codesign --verify --deep --strict --verbose=2 "$APP"
  # one submission: the signed DMG; Apple notarizes the app inside it too
  rm -f "$DMG"
  hdiutil create -volname "STAnnotator" -srcfolder "$APP" -ov -format ULMO "$DMG" \
    || hdiutil create -volname "STAnnotator" -srcfolder "$APP" -ov -format UDZO "$DMG"
  sign --sign "$DEVELOPER_ID" "$DMG"
  echo "uploading the DMG to Apple's notary service"
  xcrun notarytool submit "$DMG" --keychain-profile "$PROFILE" --output-format plist > dist/notary_submit.plist
  NOTARY_ID=$(/usr/libexec/PlistBuddy -c "Print :id" dist/notary_submit.plist)
  echo "$NOTARY_ID" > dist/notary_id.txt
fi
echo "waiting for Apple (submission $NOTARY_ID). Usually minutes; the first submissions of a new Developer ID can"
echo "take hours. Ctrl-C is safe: Apple goes on, and later"
echo "  NOTARY_ID=$NOTARY_ID bash packaging/macos/sign_and_notarize.sh"
echo "picks it up (leave dist/STAnnotator.dmg as it is)."
xcrun notarytool wait "$NOTARY_ID" --keychain-profile "$PROFILE" --timeout 2h || true
status=$(xcrun notarytool info "$NOTARY_ID" --keychain-profile "$PROFILE" | awk '/status:/ {print $2; exit}')
case "$status" in
  Accepted) ;;
  Invalid|Rejected) notary_failed "$NOTARY_ID" ;;
  *) echo "still $status after 2 h: run again later with NOTARY_ID=$NOTARY_ID"; exit 2 ;;
esac
xcrun stapler staple "$DMG"
spctl --assess --type open --context context:primary-signature -v "$DMG"
rm -f dist/notary_submit.plist dist/notary_id.txt
echo "signed and notarized: $DMG"
