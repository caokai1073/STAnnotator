# STAnnotator

Cell type annotation and analysis of spatial transcriptomics (Xenium, CosMx, MERSCOPE, Visium HD, Stereo-seq, …) on
the H&E image, with an AI agent at your side.

- **Annotate** from a marker gene list, a labeled single-cell reference, or an existing annotation.
- **Refine in plain language**: priors ("epithelial cells are only in the mucosa"), corrections, re-clustering.
  Every change is a version you can undo.
- **Analyze**: regions you draw, PCA / UMAP / Leiden, niches, DE genes, composition, cell-cell communication; one
  slice or several slices together.
- **Agent Review**: the agent reviews one part of the analysis, checks its own advice against the data, and keeps
  only what holds up.

The agent is your own [Claude Code](https://docs.claude.com/en/docs/claude-code/setup), with your own Claude account.
Your images and expression data stay on your computer (or your server).

Two ways to run it:

| | Data and computation | Window and agent |
|---|---|---|
| [A. Local computer](#a-local-computer) | your computer | your computer |
| [B. Remote server over SSH](#b-remote-server-over-ssh) | the server | your computer |

## A. Local computer

**You need**
- A Mac with Apple silicon (macOS 11+) or a PC with Windows 10 / 11 (x64), and 8 GB of RAM or more. Linux: install
  from source.
- Claude Code, installed and logged in once (`claude` works in Terminal, or in PowerShell on Windows).

**Install the app on a Mac**

Download `STAnnotator.dmg` from the [latest release](https://github.com/caokai1073/STAnnotator/releases/latest), open
it and drag **STAnnotator** to Applications. It contains everything; no Python needed.

No DMG yet? Build it on the Mac (needs [uv](https://docs.astral.sh/uv/)):

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh && source $HOME/.local/bin/env
git clone https://github.com/caokai1073/STAnnotator.git && cd STAnnotator
bash packaging/macos/build_app.sh && open dist/STAnnotator.dmg
```

**Install the app on Windows**

Download `STAnnotator-Setup.exe` from the [latest release](https://github.com/caokai1073/STAnnotator/releases/latest)
and run it. It contains everything; no Python needed.
- Choose where to install (any drive; for your user without administrator rights, or for all users) and, the first
  time, where to keep the data (samples and annotations; **Change data folder** in About moves it later).
- The installer is not signed yet: if Windows SmartScreen stops it, click **More info → Run anyway**.
- Start **STAnnotator** from the Start menu or the desktop shortcut. **Ask Agent** opens Claude Code in its own
  console window.
- Uninstalling removes the program only; the data and the settings stay.

No installer yet? Build it on the PC. Download the source (**Code → Download ZIP** on GitHub, or `git clone`), unzip
it, and in PowerShell in that folder run (nothing to install first; 5–15 minutes the first time):

```powershell
powershell -ExecutionPolicy Bypass -File packaging\windows\build.ps1
# -> dist\STAnnotator-Setup-<version>.exe, and dist\STAnnotator-windows-<version>.zip (the app folder, runs without installing)
```

**Or install from source** (Linux, developers):

```bash
git clone https://github.com/caokai1073/STAnnotator.git && cd STAnnotator
uv venv .venv --python 3.12 && uv pip install --python .venv/bin/python -e ".[app]"
.venv/bin/stagent app            # the window; --browser for a browser tab
.venv/bin/stagent make-app       # macOS: an STAnnotator launcher in ~/Applications
```

**Settings**: none. Data and annotations go to `~/stagent`, or the folder chosen in the Windows installer or with
**Change data folder** in About (`STAGENT_HOME` overrides both).

## B. Remote server over SSH

Data and computation stay on the server; your computer shows the window and runs Claude Code. Nothing is copied
to your computer except what you export.

**On the server** (Linux or macOS; no window, no root):

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh && source $HOME/.local/bin/env
git clone https://github.com/caokai1073/STAnnotator.git ~/STAnnotator && cd ~/STAnnotator
uv venv .venv --python 3.12 && uv pip install --python .venv/bin/python -e .
```

Use the same version as the app on your computer. Data goes to `~/stagent` on the server unless you choose another
folder when connecting.

**On your computer**
- The app (section A) and Claude Code, logged in.
- SSH login to the server **without a password prompt** (the app cannot type one):
  ```bash
  ssh-keygen -t ed25519            # once, if you have no key
  ssh-copy-id user@server          # then `ssh user@server` must log in directly
  ```
  On Windows (PowerShell; the OpenSSH client is part of Windows 10 / 11), instead of `ssh-copy-id`:
  ```powershell
  type $env:USERPROFILE\.ssh\id_ed25519.pub | ssh user@server "mkdir -p ~/.ssh && cat >> ~/.ssh/authorized_keys"
  ```
- Optional: a `Host` entry in `~/.ssh/config` (Windows: `%USERPROFILE%\.ssh\config`; user, port, jump host); the app
  uses it.

**Connect**: in the app click **This computer** (right panel) → **A server over SSH**, and enter
1. the server (`user@server` or a `Host` of `~/.ssh/config`),
2. the command on the server: `~/STAnnotator/.venv/bin/stagent`,
3. optionally the data folder on the server,

then **Test** and **Connect**. The app reconnects by itself at the next start; **Use this computer** switches back.

While connected, **Choose…** browses the server's files, **Export** and **Save** download to your computer, and
**Ask Agent** opens Claude Code here with its tools running on the server. No port is opened on the server: the
window reaches it through an SSH tunnel.

After updating STAnnotator on the server, reconnect once so that the server's viewer restarts with the new version.

## Getting started

1. **Import…**: Xenium (`cells.parquet` and `cell_feature_matrix.h5`; optionally `experiment.xenium` and the DAPI
   image), or h5ad / RDS (Seurat, SpatialExperiment; read without R). Optionally the H&E and its registration;
   without a registration, **Align with Agent** aligns the H&E for you.
2. **Annotate** tab: a marker list (**Markers…**: import a CSV, or **Start from an organ**: a built-in list for
   12 common organs pruned to your gene panel, which the agent can complete) or a reference, then **Annotate**.
3. Draw regions on the map; ask the agent with **Ask Agent**, e.g. "epithelial cells are only in the mucosa".
4. **Embed**, **Niches**, **Compare** and **Signaling** tabs for the analyses; **Agent Review** for a second opinion;
   **Export** for the cell table, the registration or a Jupyter notebook that reproduces the analysis.

Example data and marker lists: `examples/` (SCP2927 IBD Xenium atlas, Cadinu et al. MERFISH mouse colon).
Every control and method is described in the [user guide](docs/GUIDE.md).

## Command line

```bash
stagent import-xenium --name s1 --xenium cells.parquet cell_feature_matrix.h5 [--he he.ome.tif --registration reg.csv]
stagent import --name s1 --h5ad cells.h5ad [--he he.ome.tif --registration reg.csv --section-column donor --section 147328]
stagent viewer              # the viewer in a browser, http://localhost:8765
stagent mcp                 # the agent's tools, for Claude Code without the app
```

## Troubleshooting

- **The app does not open**: see `~/stagent/app.log` (Windows: `%LOCALAPPDATA%\STAnnotator\app.log`), or run
  `stagent app` in Terminal to see the error.
- **Ask Agent does nothing**: Claude Code must be installed and logged in (`claude` works in Terminal or PowerShell).
- **Server: Test fails**: `ssh user@server` must log in without a password, and the command path must be right
  (`ssh user@server ~/STAnnotator/.venv/bin/stagent remote-serve --check` prints JSON).
- **Few cells on the H&E**: check the registration file and that cell coordinates are in micrometers.

## Building a signed DMG (maintainers)

With an Apple Developer account, set up a Developer ID certificate and notary credentials once (see the top of
`packaging/macos/sign_and_notarize.sh`), then:

```bash
DEVELOPER_ID="Developer ID Application: Your Name (TEAMID)" bash packaging/macos/build_app.sh
```

This builds, signs and notarizes `dist/STAnnotator.dmg`, which opens on any Mac without warnings.
