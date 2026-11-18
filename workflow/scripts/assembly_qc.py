#!/usr/bin/env python3
"""QC of the filtered assembly: cells, reads, and contig lengths by sample.

Reads the FASTA written after the low-complexity filter, the candidate cell
list, and each cell's reads.fq.gz. Contig names are the SPAdes header plus
the cell id appended at assembly time (..._iN_<barcode>.<sample>).
"""

from __future__ import annotations

import argparse
import gzip
import re
import statistics
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# SPAdes writes NODE_..._i0_, _i1_, ... for isoforms. The assembler then
# appends the cell id, so the cell is whatever follows the last isoform tag.
_CELL_IN_HEADER = re.compile(r"_i\d+_(.+)$")


def cell_from_header(header: str) -> str:
    token = header[1:].split()[0]
    match = _CELL_IN_HEADER.search(token)
    if match is None:
        return ""
    return match.group(1)


def count_fastq_records(path: str) -> int:
    n_lines = 0
    with gzip.open(path, "rb") as handle:
        for _ in handle:
            n_lines += 1
    return n_lines // 4


def mean(values: list[int]) -> float:
    if not values:
        return 0.0
    return sum(values) / len(values)


def median(values: list[int]) -> float:
    if not values:
        return 0.0
    return float(statistics.median(values))


def load_cells(path: Path) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        barcode, sample = line.split("\t", 1)
        rows.append((barcode, sample, f"{barcode}.{sample}"))
    return rows


def scan_fasta(
    path: Path, known: dict[str, str]
) -> tuple[dict[str, list[int]], dict[str, int], int]:
    lengths: dict[str, list[int]] = defaultdict(list)
    contigs_by_cell: dict[str, int] = defaultdict(int)
    unknown = 0
    cell = ""
    length = 0

    def flush() -> None:
        nonlocal unknown
        if length <= 0:
            return
        sample = known.get(cell)
        if sample is None:
            unknown += 1
            return
        lengths[sample].append(length)
        contigs_by_cell[cell] += 1

    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                flush()
                cell = cell_from_header(line)
                length = 0
            else:
                length += len(line.strip())
        flush()
    return lengths, contigs_by_cell, unknown


def scan_low_complexity(path: Path, known: dict[str, str]) -> dict[str, int]:
    """Count sequences BBDuk removed, by sample. Accepts FASTQ or FASTA."""
    removed: dict[str, int] = defaultdict(int)
    if not path.is_file() or path.stat().st_size == 0:
        return removed
    with path.open() as handle:
        for line in handle:
            if not line.startswith(("@", ">")):
                continue
            sample = known.get(cell_from_header(line))
            if sample is not None:
                removed[sample] += 1
    return removed


def format_number(value: object, digits: int = 1) -> str:
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def write_report(
    samples: list[str],
    n_cells: dict[str, int],
    cells_with_contigs: dict[str, int],
    reads: dict[str, list[int]],
    missing_reads: dict[str, int],
    lengths: dict[str, list[int]],
    removed: dict[str, int],
    unknown_contigs: int,
    tsv_path: Path,
    txt_path: Path,
) -> None:
    columns = [
        "sample",
        "cells",
        "cells_with_contigs",
        "cells_without_contigs",
        "missing_reads",
        "reads",
        "mean_reads_per_cell",
        "median_reads_per_cell",
        "contigs",
        "mean_contigs_per_cell",
        "assembled_bp",
        "mean_contig_length",
        "median_contig_length",
        "contigs_gt_1000",
        "contigs_gt_1000_per_cell",
        "fraction_contigs_gt_1000",
        "low_complexity_removed",
    ]
    rows: list[dict[str, object]] = []
    for sample in samples + ["All"]:
        cell_n = n_cells[sample]
        read_counts = reads[sample]
        contig_lengths = lengths[sample]
        n_contigs = len(contig_lengths)
        assembled_bp = sum(contig_lengths)
        long_contigs = sum(1 for n in contig_lengths if n > 1000)
        with_contigs = cells_with_contigs[sample]
        rows.append(
            {
                "sample": sample,
                "cells": cell_n,
                "cells_with_contigs": with_contigs,
                "cells_without_contigs": cell_n - with_contigs,
                "missing_reads": missing_reads[sample],
                "reads": sum(read_counts),
                "mean_reads_per_cell": mean(read_counts),
                "median_reads_per_cell": median(read_counts),
                "contigs": n_contigs,
                "mean_contigs_per_cell": (n_contigs / cell_n) if cell_n else 0.0,
                "assembled_bp": assembled_bp,
                "mean_contig_length": mean(contig_lengths),
                "median_contig_length": median(contig_lengths),
                "contigs_gt_1000": long_contigs,
                "contigs_gt_1000_per_cell": (
                    (long_contigs / cell_n) if cell_n else 0.0
                ),
                "fraction_contigs_gt_1000": (
                    (long_contigs / n_contigs) if n_contigs else 0.0
                ),
                "low_complexity_removed": removed[sample],
            }
        )

    tsv_lines = ["\t".join(columns)]
    for row in rows:
        tsv_lines.append(
            "\t".join(
                format_number(
                    row[column],
                    4 if column == "fraction_contigs_gt_1000" else 1,
                )
                for column in columns
            )
        )
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    tsv_path.write_text("\n".join(tsv_lines) + "\n")

    all_row = rows[-1]
    text = [
        "Assembly QC after the low-complexity filter",
        "",
        f"Cells: {all_row['cells']}",
        (
            f"Reads per cell: mean {format_number(all_row['mean_reads_per_cell'])}, "
            f"median {format_number(all_row['median_reads_per_cell'])} "
            f"({all_row['reads']} reads)"
        ),
        (
            f"Contigs: {all_row['contigs']} "
            f"({format_number(all_row['mean_contigs_per_cell'])} per cell, "
            f"{all_row['assembled_bp']} bp)"
        ),
        (
            f"Contig length: mean {format_number(all_row['mean_contig_length'])} bp, "
            f"median {format_number(all_row['median_contig_length'])} bp"
        ),
        (
            f"Contigs > 1000 bp: {all_row['contigs_gt_1000']} "
            f"({format_number(all_row['contigs_gt_1000_per_cell'])} per cell, "
            f"{format_number(float(all_row['fraction_contigs_gt_1000']) * 100)}% of contigs)"
        ),
        f"Cells with no contigs after filtering: {all_row['cells_without_contigs']}",
        f"Low-complexity sequences removed: {all_row['low_complexity_removed']}",
    ]
    if unknown_contigs:
        text.append(f"Contigs with no candidate cell in the header: {unknown_contigs}")
    if all_row["missing_reads"]:
        text.append(f"Cells with no reads.fq.gz: {all_row['missing_reads']}")
    text.extend(["", "By sample", ""])
    text.append(
        "\t".join(
            [
                "sample",
                "cells",
                "mean_reads",
                "contigs",
                "mean_length",
                "contigs_gt_1000",
                "contigs_gt_1000_per_cell",
                "assembled_bp",
            ]
        )
    )
    for row in rows:
        text.append(
            "\t".join(
                [
                    str(row["sample"]),
                    str(row["cells"]),
                    format_number(row["mean_reads_per_cell"]),
                    str(row["contigs"]),
                    format_number(row["mean_contig_length"]),
                    str(row["contigs_gt_1000"]),
                    format_number(row["contigs_gt_1000_per_cell"]),
                    str(row["assembled_bp"]),
                ]
            )
        )
    txt_path.write_text("\n".join(text) + "\n")
    print(f"[assembly_qc] wrote {tsv_path}", flush=True)
    print(f"[assembly_qc] wrote {txt_path}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--cells-dir", required=True)
    parser.add_argument("--low-complexity", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--output-txt", required=True)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    cells = load_cells(Path(args.candidates))
    known = {cell: sample for _barcode, sample, cell in cells}
    samples: list[str] = []
    for _barcode, sample, _cell in cells:
        if sample not in samples:
            samples.append(sample)

    print(f"[assembly_qc] scanning {args.fasta}", flush=True)
    lengths, contigs_by_cell, unknown = scan_fasta(Path(args.fasta), known)
    removed = scan_low_complexity(Path(args.low_complexity), known)

    n_cells: dict[str, int] = defaultdict(int)
    cells_with_contigs: dict[str, int] = defaultdict(int)
    reads: dict[str, list[int]] = defaultdict(list)
    missing_reads: dict[str, int] = defaultdict(int)
    jobs: list[tuple[str, str]] = []
    cells_dir = Path(args.cells_dir)
    for _barcode, sample, cell in cells:
        n_cells[sample] += 1
        n_cells["All"] += 1
        if contigs_by_cell.get(cell, 0) > 0:
            cells_with_contigs[sample] += 1
            cells_with_contigs["All"] += 1
        reads_path = cells_dir / cell / "reads.fq.gz"
        if not reads_path.is_file() or reads_path.stat().st_size == 0:
            missing_reads[sample] += 1
            missing_reads["All"] += 1
            reads[sample].append(0)
            reads["All"].append(0)
            continue
        jobs.append((sample, str(reads_path)))

    print(f"[assembly_qc] counting reads for {len(jobs)} cells", flush=True)
    done = 0
    workers = max(1, args.threads)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(count_fastq_records, path): sample for sample, path in jobs
        }
        for future in as_completed(futures):
            sample = futures[future]
            n_reads = future.result()
            reads[sample].append(n_reads)
            reads["All"].append(n_reads)
            done += 1
            if done % 2000 == 0:
                print(f"[assembly_qc] reads counted {done}/{len(jobs)}", flush=True)

    all_lengths: list[int] = []
    all_removed = 0
    for sample in samples:
        all_lengths.extend(lengths[sample])
        all_removed += removed[sample]
    lengths["All"] = all_lengths
    removed["All"] = all_removed

    write_report(
        samples,
        n_cells,
        cells_with_contigs,
        reads,
        missing_reads,
        lengths,
        removed,
        unknown,
        Path(args.output_tsv),
        Path(args.output_txt),
    )


if __name__ == "__main__":
    main()
