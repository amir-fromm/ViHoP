#!/usr/bin/env python3
"""Fail unless every split has a DIAMOND virus summary."""

from __future__ import annotations

import sys
from pathlib import Path


def main() -> None:
    n = int(sys.argv[1])
    if n < 1:
        raise SystemExit(f"batch count must be >= 1, got {n}")
    parts_dir = Path("results/virus/parts")
    missing: list[str] = []
    for i in range(n):
        path = parts_dir / f"batch_{i:04d}.tsv"
        if not path.is_file():
            missing.append(str(path))
    if missing:
        print(
            f"ERROR: {len(missing)} virus batch summaries missing. "
            f"All {n} splits must finish DIAMOND before the tables are merged.",
            file=sys.stderr,
        )
        for path in missing:
            print(path, file=sys.stderr)
        sys.exit(1)
    extra = sorted(
        p.name
        for p in parts_dir.glob("batch_*.tsv")
        if p.stem not in {f"batch_{i:04d}" for i in range(n)}
    )
    if extra:
        print(
            f"ERROR: {len(extra)} virus part(s) on disk are outside 0..{n - 1}. "
            f"params.host.split_fasta/barrnap_batches={n} is too small "
            f"(e.g. {extra[0]}). Fix the config before merging.",
            file=sys.stderr,
        )
        sys.exit(1)
    print(f"[verify] all {n} splits finished the virus DIAMOND search", flush=True)


if __name__ == "__main__":
    main()
