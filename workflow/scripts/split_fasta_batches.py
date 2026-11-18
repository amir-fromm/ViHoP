#!/usr/bin/env python3
"""Round-robin split of a FASTA into N batch files."""

from __future__ import annotations

import sys
from pathlib import Path


def main() -> None:
    src, outdir_s, n_s = sys.argv[1:4]
    n = int(n_s)
    if n < 1:
        raise SystemExit(f"batch count must be >= 1, got {n}")
    outdir = Path(outdir_s)
    outdir.mkdir(parents=True, exist_ok=True)
    handles = [(outdir / f"batch_{i:04d}.fasta").open("w") for i in range(n)]
    index = 0
    seen = 0
    header = None
    body: list[str] = []

    def flush() -> None:
        nonlocal index, seen, header, body
        if header is None:
            return
        handle = handles[index % n]
        handle.write(header if header.endswith("\n") else header + "\n")
        handle.writelines(body)
        if body and not body[-1].endswith("\n"):
            handle.write("\n")
        index += 1
        seen += 1
        header = None
        body = []

    with open(src) as handle:
        for line in handle:
            if line.startswith(">"):
                flush()
                header = line
            elif header is not None:
                body.append(line)
        flush()
    for handle in handles:
        handle.close()
    print(f"[split] sequences={seen} batches={n} dir={outdir}", flush=True)


if __name__ == "__main__":
    main()
