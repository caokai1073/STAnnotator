# Build the standalone STAnnotator app on Windows 10 / 11 (x64): dist\STAnnotator\STAnnotator.exe (the window) and
# stagent.exe beside it (the console program Claude Code starts), zipped as dist\STAnnotator-windows-<version>.zip.
#
# From the unzipped STAnnotator folder, in PowerShell:
#     powershell -ExecutionPolicy Bypass -File packaging\windows\build.ps1
#
# It needs Python 3.12: uv (https://docs.astral.sh/uv/) when it is installed, else the Python launcher (py -3.12)
# of a python.org installation. The build goes into its own environment (.build-venv), without PyTorch, and takes
# 5 - 15 minutes the first time (downloads about 1 GB of packages).
$ErrorActionPreference = "Stop"
Set-Location (Join-Path $PSScriptRoot "..\..")
$root = Get-Location
$venv = Join-Path $root ".build-venv"
$py = Join-Path $venv "Scripts\python.exe"

if (Get-Command uv -ErrorAction SilentlyContinue) {
    Write-Host "== Python 3.12 environment (uv)"
    uv venv --allow-existing $venv --python 3.12
    uv pip install --python $py --reinstall-package stagent ".[app]" pyinstaller
} elseif (Get-Command py -ErrorAction SilentlyContinue) {
    Write-Host "== Python 3.12 environment (py -3.12)"
    if (-not (Test-Path $py)) { py -3.12 -m venv $venv }
    & $py -m pip install --upgrade pip
    & $py -m pip install --force-reinstall --no-deps .
    & $py -m pip install ".[app]" pyinstaller
} else {
    Write-Error ("Python 3.12 was not found. Install uv (https://docs.astral.sh/uv/getting-started/installation/) " +
                 "or Python 3.12 from https://www.python.org/downloads/windows/ (tick 'Add python.exe to PATH'), " +
                 "then run this script again.")
}

Write-Host "== PyInstaller"
& $py -m PyInstaller --noconfirm --clean packaging\stannotator.spec
if ($LASTEXITCODE -ne 0) { throw "PyInstaller failed" }

$app = Join-Path $root "dist\STAnnotator"
Write-Host "== check: stagent.exe answers"
& (Join-Path $app "stagent.exe") --help | Select-Object -First 3
if ($LASTEXITCODE -ne 0) { throw "the built stagent.exe does not start" }

$version = (Select-String -Path stagent\__init__.py -Pattern '__version__ = "([^"]+)"').Matches[0].Groups[1].Value
$zip = Join-Path $root "dist\STAnnotator-windows-$version.zip"
if (Test-Path $zip) { Remove-Item $zip }
Write-Host "== zip"
Compress-Archive -Path $app -DestinationPath $zip -CompressionLevel Optimal
$size = "{0:N0} MB" -f ((Get-ChildItem $app -Recurse | Measure-Object Length -Sum).Sum / 1MB)
Write-Host ""
Write-Host "built: $(Get-Content stagent\data\build.txt -ErrorAction SilentlyContinue)  version $version"
Write-Host "app:   $app\STAnnotator.exe  ($size)"
Write-Host "zip:   $zip"
Write-Host "Start it with a double-click on STAnnotator.exe (Windows may first ask: More info -> Run anyway)."
