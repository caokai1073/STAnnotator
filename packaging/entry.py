"""Entry point of the standalone app: no arguments start the desktop app; otherwise the celltypeagent subcommands
(mcp for Claude, viewer, import, ...)."""
import multiprocessing
import sys

if __name__ == "__main__":
    multiprocessing.freeze_support()
    args = [a for a in sys.argv[1:] if not a.startswith("-psn_")]  # Finder may pass a process serial number
    from celltypeagent.cli import main
    main(args or ["app"])
