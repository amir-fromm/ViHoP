#!/usr/bin/env python3
"""Fail unless every batch has Barrnap outputs and a PR2 BLAST table."""

from __future__ import annotations

import sys
from pathlib import Path


def main() -> None:
    n = int(sys.argv[1])
    if n < 1:
        raise SystemExit(f"batch count must be >= 1, got {n}")
    missing: list[str] = []
    for i in range(n):
        batch = f"batch_{i:04d}"
        needed = (
            Path(f"results/host/barrnap_parts/{batch}.rrna.fa"),
            Path(f"results/host/barrnap_parts/{batch}.rrna.gff"),
            Path(f"results/host/barrnap_parts/{batch}.18s.fasta"),
            Path(f"results/host/blast_parts/{batch}.PR2.tsv"),
        )
        for path in needed:
            if not path.is_file():
                missing.append(str(path))
    if missing:
        print(
            f"ERROR: {len(missing)} files missing. "
            f"All {n} batches must finish Barrnap and PR2 BLAST "
            "before summary.",
            file=sys.stderr,
        )
        for path in missing:
            print(path, file=sys.stderr)
        sys.exit(1)
    print(
        f"[verify] all {n} batches finished Barrnap and PR2 BLAST",
        flush=True,
    )


if __name__ == "__main__":
    main()
