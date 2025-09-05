"""Ad-hoc test for Pantheon-CLI access prompt flow (no network).

Runs the REPL UI tool result printer with a fake Domain Research result
containing `access_requests`. Verifies auto-open behavior and interactive
prompt wiring without requiring full REPL or network access.

Usage:
  python pantheon-cli/scripts/test_cli_access_flow.py --mode oa
  python pantheon-cli/scripts/test_cli_access_flow.py --mode proxy
  python pantheon-cli/scripts/test_cli_access_flow.py --mode none

Env:
  ALLOW_BROWSER_OPEN=1 to enable opening links (monkeypatched in test).
"""

from __future__ import annotations

import argparse
import io
import os
import sys
import types
from typing import List


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["oa", "proxy", "none"], default="none")
    args = parser.parse_args()

    # Import UI
    # Prefer local package path
    sys.path.insert(0, "pantheon-cli")
    from pantheon_cli.repl.ui import ReplUI

    # Prepare a fake result (as returned by domain_research.run_research)
    result = {
        "report": "# Demo Report\n...",
        "backend_used": "web:pubmed",
        "used_demo": False,
        "access_requests": [
            {
                "title": "Open-access item",
                "url": "https://journals.example.org/item",
                "open_access_url": "https://oa.example.org/item.pdf",
                "proxy_url": "https://libproxy.example.edu/login?url=https%3A%2F%2Fjournals.example.org%2Fitem",
            },
            {
                "title": "Proxy-only item",
                "url": "https://publisher.example.com/locked",
                "proxy_url": "https://libproxy.example.edu/login?url=https%3A%2F%2Fpublisher.example.com%2Flocked",
                "access_required": True,
            },
            {
                "title": "OA-only item",
                "url": "https://preprint.example.net/article",
                "open_access_url": "https://preprint.example.net/article.pdf",
            },
            {
                "title": "No links item (should be skipped)",
                "url": "https://unknown.example.net/",
            },
        ],
    }

    # Prepare UI with captured console
    ui = ReplUI()
    out = io.StringIO()
    from rich.console import Console

    ui.console = Console(file=out, force_terminal=False, no_color=True)

    # Monkeypatch webbrowser.open to capture calls
    import webbrowser as _wb
    # Monkeypatch Prompt.ask to avoid interactive input
    from rich import prompt as _rich_prompt

    opened: List[str] = []

    def fake_open(url: str, *a, **k):
        opened.append(url)
        return True

    # Set envs
    os.environ["ALLOW_BROWSER_OPEN"] = "1"
    os.environ["ACCESS_OPEN_MODE"] = args.mode

    # Apply monkeypatch via assignment
    _wb_open_orig = _wb.open
    _wb.open = fake_open  # type: ignore
    _orig_ask = _rich_prompt.Prompt.ask
    def _fake_ask(*a, **k):
        return 's'
    _rich_prompt.Prompt.ask = _fake_ask  # type: ignore

    try:
        ui.print_tool_result("domain_research.run_research", result)
        # Call prompt explicitly to ensure coverage
        ui._prompt_access_requests(result["access_requests"])  # type: ignore
    finally:
        _rich_prompt.Prompt.ask = _orig_ask  # restore
        _wb.open = _wb_open_orig  # restore

    # Emit summary
    print("=== Captured Output ===")
    print(out.getvalue())
    print("=== Opened URLs (mode=%s) ===" % args.mode)
    for u in opened:
        print(u)
    print("COUNT:", len(opened))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
