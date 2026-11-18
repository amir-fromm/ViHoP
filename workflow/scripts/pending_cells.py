#!/usr/bin/env python3
"""Write or delete the pending extraction and assembly cell lists.

    pending_cells.py write  WORKDIR
    pending_cells.py delete WORKDIR

write scans WORKDIR/results/cells (if it exists) against candidate_cells.tsv.
Cells with a non-empty reads.fq.gz are left out of results/tmp_extract.tsv.
Cells with non-empty transcripts.fasta and transcripts.edit.fasta are left
out of results/tmp_assemble.tsv. delete removes those lists when the run ends.
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import delete_pending_cell_lists, write_pending_cell_lists


def write_lists(workdir: Path) -> None:
    results = workdir / "results"
    candidates = results / "candidate_cells.tsv"
    if not candidates.is_file():
        print(f"[pending] no {candidates}", file=sys.stderr)
        sys.exit(1)
    extract_rows, assemble_rows = write_pending_cell_lists(
        candidates,
        results / "cells",
        results / "tmp_extract.tsv",
        results / "tmp_assemble.tsv",
        results / "tmp_pending.ready",
    )
    extract_by = Counter(sample for _, sample in extract_rows)
    assemble_by = Counter(sample for _, sample in assemble_rows)
    print(
        f"[pending] tmp_extract={len(extract_rows)} "
        f"tmp_assemble={len(assemble_rows)}",
        flush=True,
    )
    for sample in sorted(set(extract_by) | set(assemble_by)):
        print(
            f"[pending] {sample} extract={extract_by[sample]} "
            f"assemble={assemble_by[sample]}",
            flush=True,
        )


def main() -> None:
    if len(sys.argv) != 3 or sys.argv[1] not in {"write", "delete"}:
        print(f"Usage: {sys.argv[0]} write|delete WORKDIR", file=sys.stderr)
        sys.exit(2)
    workdir = Path(sys.argv[2])
    if sys.argv[1] == "write":
        write_lists(workdir)
        return
    delete_pending_cell_lists(workdir / "results")
    print("[pending] removed tmp_extract and tmp_assemble", flush=True)


if __name__ == "__main__":
    main()
