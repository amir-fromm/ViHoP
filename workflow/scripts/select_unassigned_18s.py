#!/usr/bin/env python3
"""Write 18S sequences that have no PR2 hit at the identity and length cutoffs."""

from __future__ import annotations

import sys
from pathlib import Path


def assigned_ids(tsv: Path, identity: float, min_length: float) -> set[str]:
    """Queries whose best PR2 hit (by bit score) passes both cutoffs."""
    best: dict[str, tuple[float, float, float]] = {}
    if not tsv.is_file() or tsv.stat().st_size == 0:
        return set()
    with tsv.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            try:
                pident = float(fields[2])
                length = float(fields[3])
                bitscore = float(fields[11])
            except ValueError:
                continue
            current = best.get(fields[0])
            if current is None or bitscore > current[0]:
                best[fields[0]] = (bitscore, pident, length)
    return {
        query
        for query, (_score, pident, length) in best.items()
        if pident >= identity and length >= min_length
    }


def main() -> None:
    fasta, tsv, out, identity_s, length_s = sys.argv[1:6]
    identity = float(identity_s)
    min_length = float(length_s)
    keep_out = assigned_ids(Path(tsv), identity, min_length)
    written = 0
    total = 0
    header = None
    body: list[str] = []
    dest = Path(out)
    dest.parent.mkdir(parents=True, exist_ok=True)

    def flush(handle) -> None:
        nonlocal written, header, body
        if header is None:
            return
        query = header[1:].split()[0]
        if query not in keep_out:
            handle.write(header if header.endswith("\n") else header + "\n")
            handle.writelines(body)
            if body and not body[-1].endswith("\n"):
                handle.write("\n")
            written += 1
        header = None
        body = []

    with dest.open("w") as handle, open(fasta) as source:
        for line in source:
            if line.startswith(">"):
                flush(handle)
                header = line
                total += 1
            elif header is not None:
                body.append(line)
        flush(handle)
    print(
        f"[metapr2-query] sequences={total} assigned_by_pr2={total - written} "
        f"for_metapr2={written} identity>={identity} length>={min_length}",
        flush=True,
    )


if __name__ == "__main__":
    main()
