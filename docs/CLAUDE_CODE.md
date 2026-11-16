# Setting up Claude Code for STAnnotator

STAnnotator's agent is **Claude Code**, Anthropic's command-line assistant, running on your own computer with your
own Claude account. The app computes every result itself; Claude Code reads those results through the app's tools and
proposes changes, which wait in the app until you accept them.

You need Claude Code for the agent features only:

- the **Agent** button (talking to Claude about your sample),
- **Refine cell types** with **Agent review**,
- **Review** of an annotation, of the registration or of an analysis,
- aligning an H&E image without a registration file (**Align with Agent**).

Without Claude Code, or offline, everything else works: import, registration with a file, marker-based annotation,
regions, priors you write yourself, the analyses and the exports.

Setting it up takes about 10 minutes and is done once per computer.

## 1. A Claude plan that includes Claude Code

Claude Code needs one of these accounts:

- a paid Claude plan: **Pro**, **Max**, **Ultra**, **Team** or **Enterprise** (see
  [claude.com/pricing](https://claude.com/pricing)), or
- a **Claude Console** account (pay per use with an API key, at [platform.claude.com](https://platform.claude.com)).

The free Claude plan does not include Claude Code. Claude is available only in
[Anthropic's supported countries](https://www.anthropic.com/supported-countries).

How much you can use depends on the plan. A **Refine cell types** round or a **Review** is a few minutes of the
agent's work; Pro is enough to try STAnnotator, and Max gives more room for many rounds a day. When a limit is
reached, Claude Code says so in its window and when it resets.

## 2. Install Claude Code

Claude Code runs on macOS 13 or later, Windows 10 (1809) or later and common Linux distributions, with 4 GB of memory
or more and an internet connection.

### Mac

1. Open **Terminal**: press <kbd>⌘</kbd> + <kbd>Space</kbd>, type `Terminal`, press <kbd>Return</kbd>.
2. Paste this line and press <kbd>Return</kbd>:

   ```bash
   curl -fsSL https://claude.ai/install.sh | bash
   ```

   It shows no progress while it downloads; wait until the prompt comes back.
3. Close Terminal and open a new Terminal window, then check:

   ```bash
   claude --version
   ```

   It prints a version number, such as `2.1.211 (Claude Code)`.

With [Homebrew](https://brew.sh) you can instead run `brew install --cask claude-code` (Homebrew installs do not
update themselves: run `brew upgrade claude-code` from time to time).

### Windows

1. Open **PowerShell**: click Start, type `PowerShell`, press <kbd>Enter</kbd>. The prompt starts with `PS C:\`.
   No administrator rights are needed.
2. Paste this line and press <kbd>Enter</kbd>:

   ```powershell
   irm https://claude.ai/install.ps1 | iex
   ```

   In the older **Command Prompt** (a prompt without `PS`), use this line instead:

   ```bat
   curl -fsSL https://claude.ai/install.cmd -o install.cmd && install.cmd && del install.cmd
   ```

3. Close the window, open a new PowerShell window, then check:

   ```powershell
   claude --version
   ```

Or install it with WinGet: `winget install Anthropic.ClaudeCode` (update it with
`winget upgrade Anthropic.ClaudeCode`). [Git for Windows](https://git-scm.com/downloads/win) is optional: Claude Code
uses PowerShell without it.

### Linux server

The same line as on a Mac, in a terminal on the server:

```bash
curl -fsSL https://claude.ai/install.sh | bash
```

You need Claude Code on a server only for the server's web viewer and its chat box. When you work from the app on
your computer with **Local → A server over SSH**, Claude Code runs on your computer (its tools run on the server
through `ssh`): install it on your computer, not on the server.

### If `claude` is not found

If the new terminal says `command not found: claude` (Mac) or `'claude' is not recognized` (Windows), the folder
it was installed in is not on your PATH yet. Open one more new terminal window; if it is still not found, follow
[Fix your PATH](https://code.claude.com/docs/en/troubleshoot-install.md#command-not-found-claude-after-installation).
`claude doctor` checks the installation and says what is wrong.

The installer keeps Claude Code up to date by itself.

## 3. Log in once

In the terminal, run:

```bash
claude
```

The first time, Claude Code asks for a color theme, then opens your browser to log in. Log in with the account of
your plan (or choose the Console account), allow the access, and go back to the terminal. When Claude Code shows its
prompt, it is ready: type `/exit` (or press <kbd>Ctrl</kbd> + <kbd>C</kbd> twice) to leave it.

You stay logged in; you do not need to do this again.

## 4. Connect it to STAnnotator

There is nothing to configure: the app finds Claude Code and connects it to its tools by itself.

1. Start STAnnotator. Under **?** → **About**, the line **Claude Code** says **found** with where it is. If it says
   **not found**, quit the app and start it again (it looks for Claude Code when it starts).
2. Open a sample and click **Agent** (on the bar at the top of the map).
   - **Mac**: the first time, macOS asks whether STAnnotator may control **Terminal**. Click **OK** (or **Allow**). (If you clicked
     **Don't Allow**, turn it on in **System Settings → Privacy & Security → Automation → STAnnotator → Terminal**.)
   - **Windows**: Claude Code opens in a console window of its own.
3. The first time, Claude Code asks whether you trust the folder it opens in (`claude_workspace` in the data folder,
   e.g. `~/stagent/claude_workspace`). Choose **Yes**: the folder holds only the connection to the app.
4. Claude Code starts connected to the app and to the sample you are looking at. Write to it in your own words, for
   example "Review my annotation" or "Epithelial cells are only in the mucosa". The app's tools are allowed in that
   folder; for anything else Claude Code asks you first.

**Refine cell types** (after a marker annotation) sends its request to Claude Code by itself when you
click **Agent review**; its proposal appears in the Refine window.

The box next to **Agent** chooses the model (**Fable**, **Opus**, **Sonnet**, **Haiku**, or **Default**, Claude
Code's own setting). In a window already open, type `/model sonnet` and so on.

## 5. Common problems

| What you see | What to do |
|---|---|
| **About** says Claude Code **not found** | In a new terminal, `claude --version` must work. Then quit and restart STAnnotator. |
| Clicking **Agent** opens nothing (Mac) | Allow STAnnotator to control Terminal: **System Settings → Privacy & Security → Automation**. |
| Claude Code asks you to log in every time, or says the account has no access | Log in with a paid plan or a Console account (step 1 and 3). The free plan does not include Claude Code. |
| "No internet" in the app, or Claude Code cannot connect | Claude Code needs the internet. Behind a company proxy, see [network configuration](https://code.claude.com/docs/en/network-config.md). |
| A usage limit message in the Claude Code window | Wait until the time it gives, or use a plan with more usage. Nothing is lost: click **Agent review** again afterwards. |
| **Refine cell types** waits more than 10 minutes | Look at the Claude Code window: it shows the error there (no internet, not logged in, a limit reached) or is still working. |

To open Claude Code connected to the app by hand (Mac or Linux):

```bash
cd ~/stagent/claude_workspace && claude --mcp-config stagent_mcp.json
```

(Use your data folder if it is not `~/stagent`; on Windows, `cd /d` into it in Command Prompt.)

## What Claude sees

Claude Code works under your own Claude account and its terms. It receives what you write and what the app's tools
return to it: numbers, tables, the marker list, and figures of the tissue it asks for. The app does not upload your
data files; they stay in your data folder.

More: [Claude Code setup](https://code.claude.com/docs/en/setup) and
[troubleshooting](https://code.claude.com/docs/en/troubleshoot-install.md) in Anthropic's documentation, and the
[STAnnotator user guide](GUIDE.md).
