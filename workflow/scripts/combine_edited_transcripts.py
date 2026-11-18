#!/usr/bin/env python3
"""Concatenate per-cell transcripts.edit.fasta listed in candidate_cells.tsv."""

from __future__ import annotations

import sys
from pathlib import Path


def main() -> None:
    candidates, output = sys.argv[1:3]
    dest = Path(output)
    dest.parent.mkdir(parents=True, exist_ok=True)
    missing: list[str] = []
    written = 0
    with dest.open("w") as handle:
        for line in open(candidates):
            if not line.strip():
                continue
            barcode, sample = line.rstrip("\n").split("\t")
            src = Path("results/cells") / f"{barcode}.{sample}" / "transcripts.edit.fasta"
            if not src.is_file() or src.stat().st_size == 0:
                missing.append(str(src))
                continue
            handle.write(src.read_text())
            written += 1
    if missing:
        print(
            f"ERROR: {len(missing)} edited transcripts missing (first: {missing[0]})",
            file=sys.stderr,
        )
        sys.exit(1)
    print(f"[combine] cells={written} -> {dest}", flush=True)


if __name__ == "__main__":
    main()
