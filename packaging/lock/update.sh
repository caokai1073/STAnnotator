#!/bin/bash
# The release builds install exactly the packages of packaging/lock/{mac,win}.txt (versions and hashes), so a build of
# the same code gives the same app whenever it runs. A new version of a package, released the day of a build, broke
# 0.2.55 that way (imagecodecs 2026.11.11: the app could no longer read or write JPEG tiles).
#
# Updating the packages is a step of its own: run this (on any computer with uv), build and test the app (the builds
# run "stagent selfcheck" on it; run the smoke tests from source with the new versions), and commit the two files.
# With a date (e.g. 2026-12-01T00:00:00Z) it takes the packages released before then; without one, the newest.
# pyproject.toml keeps the ranges (pins where versions change results); the lock files are the exact versions in them.
set -euo pipefail
cd "$(dirname "$0")/../.."
EXCLUDE=()
[ -n "${1:-}" ] && EXCLUDE=(--exclude-newer "$1")
for p in mac:aarch64-apple-darwin win:x86_64-pc-windows-msvc; do
  uv pip compile pyproject.toml packaging/lock/build.in --extra app --python-version 3.12 \
    --python-platform "${p#*:}" "${EXCLUDE[@]}" --generate-hashes --no-header -q -o "packaging/lock/${p%%:*}.txt"
done
# (MLX is not in mac.txt: its marker needs the Mac's kernel version; build_app.sh installs the version pyproject pins)
git diff --stat -- packaging/lock/
