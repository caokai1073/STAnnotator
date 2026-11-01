# Build the standalone STAnnotator app on Windows 10 / 11 (x64): dist\STAnnotator\STAnnotator.exe (the window) and
# stagent.exe beside it (the console program Claude Code starts), zipped as dist\STAnnotator-windows-<version>.zip.
#
# From the unzipped STAnnotator folder, in PowerShell:
#     powershell -ExecutionPolicy Bypass -File packaging\windows\build.ps1
#
# Nothing needs to be installed first: it uses uv (https://docs.astral.sh/uv/) when it is installed, else downloads uv
# into .build-tools\ of this folder (checked against its published SHA-256; nothing is installed on the system), and
# uv downloads Python 3.12 for the build. The build goes into its own environment (.build-venv), without PyTorch, and
# takes 5 - 15 minutes the first time (downloads about 1 GB). The app it makes contains Python: whoever runs
# STAnnotator.exe needs nothing installed.
$ErrorActionPreference = "Stop"
Set-Location (Join-Path $PSScriptRoot "..\..")
$root = Get-Location
$venv = Join-Path $root ".build-venv"
$py = Join-Path $venv "Scripts\python.exe"

$uv = (Get-Command uv -ErrorAction SilentlyContinue).Source
if (-not $uv) {
    # uv for this build only, in .build-tools\ (no installation, no change to PATH)
    $tools = Join-Path $root ".build-tools"
    $uv = Join-Path $tools "uv.exe"
    if (-not (Test-Path $uv)) {
        Write-Host "== downloading uv (Python and packages for the build come through it)"
        New-Item -ItemType Directory -Force $tools | Out-Null
        $url = "https://github.com/astral-sh/uv/releases/latest/download/uv-x86_64-pc-windows-msvc.zip"
        $zipUv = Join-Path $tools "uv.zip"
        [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
        Invoke-WebRequest -Uri $url -OutFile $zipUv -UseBasicParsing
        $shaFile = Join-Path $tools "uv.zip.sha256"
        Invoke-WebRequest -Uri "$url.sha256" -OutFile $shaFile -UseBasicParsing   # (as a file: it may arrive as bytes)
        $want = ((Get-Content $shaFile -Raw).Trim() -split "\s+")[0].ToLower()
        Remove-Item $shaFile
        $have = (Get-FileHash $zipUv -Algorithm SHA256).Hash.ToLower()
        if ($want -ne $have) { Remove-Item $zipUv; throw "uv download: checksum mismatch ($have, expected $want)" }
        Expand-Archive -Path $zipUv -DestinationPath $tools -Force
        Remove-Item $zipUv
        if (-not (Test-Path $uv)) { $uv = (Get-ChildItem $tools -Recurse -Filter uv.exe | Select-Object -First 1).FullName }
    }
}
Write-Host "== Python 3.12 environment (uv: $uv)"
& $uv venv --allow-existing $venv --python 3.12      # uv downloads Python 3.12 when the computer has none
if ($LASTEXITCODE -ne 0) { throw "uv could not create the Python 3.12 environment" }
& $uv pip install --python $py --reinstall-package stagent ".[app]" pyinstaller
if ($LASTEXITCODE -ne 0) { throw "installing the packages failed" }

Write-Host "== PyInstaller"
& $py -m PyInstaller --noconfirm --clean packaging\stannotator.spec
if ($LASTEXITCODE -ne 0) { throw "PyInstaller failed" }

$app = Join-Path $root "dist\STAnnotator"
Write-Host "== check: stagent.exe answers"
# the whole output first: Select-Object -First in the pipe would stop the program early and spoil its exit code
$out = & (Join-Path $app "stagent.exe") --help 2>&1
$code = $LASTEXITCODE
$out | Select-Object -First 3
if ($code -ne 0) { throw "the built stagent.exe does not start (exit code $code)" }

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
