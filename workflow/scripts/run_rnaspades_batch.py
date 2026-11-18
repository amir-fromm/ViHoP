#!/usr/bin/env python3
"""Assemble many cells in one job. Skip cells that already have transcripts."""

from __future__ import annotations

import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def edit_headers(transcripts: Path, edited: Path, cell: str) -> None:
    lines: list[str] = []
    for line in transcripts.read_text().splitlines():
        if line.startswith(">"):
            lines.append(f"{line}_{cell}")
        else:
            lines.append(line)
    text = "\n".join(lines)
    if text:
        text += "\n"
    edited.write_text(text)


def assemble_one(barcode: str, sample: str, cell_script: str) -> str | None:
    cell = f"{barcode}.{sample}"
    folder = Path("results/cells") / cell
    reads = folder / "reads.fq.gz"
    transcripts = folder / "transcripts.fasta"
    edited = folder / "transcripts.edit.fasta"
    reads_newer = (
        transcripts.is_file()
        and reads.is_file()
        and reads.stat().st_mtime > transcripts.stat().st_mtime
    )
    if (
        not transcripts.is_file()
        or transcripts.stat().st_size == 0
        or reads_newer
    ):
        if not reads.is_file() or reads.stat().st_size == 0:
            return f"missing reads: {reads}"
        subprocess.run(
            ["bash", cell_script, str(reads), str(transcripts), "1"],
            check=True,
        )
    edit_headers(transcripts, edited, cell)
    return None


def main() -> None:
    manifest, doneflag, threads_s, cell_script = sys.argv[1:5]
    threads = max(1, int(threads_s))
    rows: list[tuple[str, str]] = []
    for line in open(manifest):
        if not line.strip():
            continue
        barcode, sample = line.rstrip("\n").split("\t")
        rows.append((barcode, sample))

    errors: list[str] = []
    with ThreadPoolExecutor(max_workers=threads) as pool:
        futures = [
            pool.submit(assemble_one, barcode, sample, cell_script)
            for barcode, sample in rows
        ]
        for future in as_completed(futures):
            try:
                err = future.result()
            except subprocess.CalledProcessError as exc:
                errors.append(f"spades failed ({exc.returncode}): {exc.cmd}")
            else:
                if err:
                    errors.append(err)

    if errors:
        print(f"[assemble-batch] {len(errors)} failures", file=sys.stderr)
        for err in errors[:30]:
            print(err, file=sys.stderr)
        sys.exit(1)

    done = Path(doneflag)
    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text(f"cells={len(rows)}\n")
    print(f"[assemble-batch] done cells={len(rows)} -> {done}", flush=True)


if __name__ == "__main__":
    main()
