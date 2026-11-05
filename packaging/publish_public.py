"""Fill the public repository (github.com/caokai1073/STAnnotator: website, README, user guide, releases; no source
code) from this one. The code lives in the private repository; the public one is a separate checkout, by default
./public (ignored here).

    python packaging/publish_public.py [public_checkout]

Copies the website (docs/index.html, docs/.nojekyll, the assets it uses), packaging/public/README.md as README.md,
LICENSE, and docs/GUIDE.md without what only works from the source (command-line imports, the server installation,
other clients, the code layout, the bundled examples). Then commit and push in the public checkout.
"""
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ASSETS = ["logo.png", "favicon.png", "apple-touch-icon.png", "social-preview.png", "app-overview.jpg",
          "mucosa-detail.jpg", "overview.png"]
SOURCE_ONLY = re.compile(r"git clone|\.venv/|uv pip|uv venv|python -m stagent|examples/")


def drop_section(text, heading):
    """Remove a section from its heading to the next heading of the same or a higher level."""
    level = heading.split(" ")[0]
    start = text.find(heading + "\n")
    if start < 0:
        return text
    rest = text[start + len(heading) + 1:]
    m = re.search(r"^#{1,%d} " % len(level), rest, flags=re.M)
    return text[:start] + (rest[m.start():] if m else "")


def public_guide(text):
    for heading in ("### Example: the IBD Xenium atlas (SCP2927)", "## Other clients", "## Layout"):
        text = drop_section(text, heading)
    # command-line imports
    text = re.sub(r"From the command line \(`--he` and `--registration` are optional\):\n\n```bash\n.*?```\n\n",
                  "", text, flags=re.S)
    # the server installation: available on request
    text = re.sub(r"On the server, install STAnnotator once[^\n]*\n\n```bash\n.*?\n```\n",
                  "The server side is available on request: "
                  "[open an issue](https://github.com/caokai1073/STAnnotator/issues).\n", text, flags=re.S)
    text = text.replace(" (or the web viewer: `.venv/bin/stagent viewer`, then http://localhost:8765)", "")
    text = text.replace(" by `python -m stagent.lr_resource`;", ";")
    text = text.replace("; run `.venv/bin/stagent app` in a terminal to see\n  errors, or `--browser` to use the browser",
                        "")
    return text


def main():
    dest = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "public"
    (dest / "docs" / "assets").mkdir(parents=True, exist_ok=True)
    shutil.copy(ROOT / "packaging" / "public" / "README.md", dest / "README.md")
    shutil.copy(ROOT / "LICENSE", dest / "LICENSE")
    shutil.copy(ROOT / "docs" / "index.html", dest / "docs" / "index.html")
    (dest / "docs" / ".nojekyll").touch()
    for name in ASSETS:
        shutil.copy(ROOT / "docs" / "assets" / name, dest / "docs" / "assets" / name)
    guide = public_guide((ROOT / "docs" / "GUIDE.md").read_text(encoding="utf-8"))
    (dest / "docs" / "GUIDE.md").write_text(guide, encoding="utf-8")
    left = [f"{p.relative_to(dest)}:{i}: {line.strip()[:90]}"
            for p in [dest / "README.md", dest / "docs" / "GUIDE.md", dest / "docs" / "index.html"]
            for i, line in enumerate(p.read_text(encoding="utf-8").splitlines(), 1) if SOURCE_ONLY.search(line)]
    print(f"public repository filled: {dest}")
    if left:
        print("still mentions the source:\n  " + "\n  ".join(left))


if __name__ == "__main__":
    main()
