#!/usr/bin/env python3
"""Record one finished step in results/finished/steps.yaml."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import mark_step_finished  # noqa: E402


def main() -> None:
    if len(sys.argv) != 2 or not sys.argv[1].strip():
        raise SystemExit("usage: mark_step_finished.py STEP_NAME")
    mark_step_finished(sys.argv[1].strip())
    print(f"[finished] {sys.argv[1].strip()}", flush=True)


if __name__ == "__main__":
    main()
