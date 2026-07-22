#!/usr/bin/env python3
"""Extract reads for candidate cells from a Cell Ranger BAM (one pass)."""

import argparse
import gzip
import os
import sys
import time

import pysam

parser = argparse.ArgumentParser()
parser.add_argument("--bam", required=True)
parser.add_argument("--cells", required=True)
parser.add_argument("--outdir", required=True)
parser.add_argument("--threads", type=int, default=2)
parser.add_argument(
    "--progress-every",
    type=int,
    default=5_000_000,
    help="Print progress every N BAM records (0 disables).",
)
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

writers = {}
written = {}

with open(args.cells) as handle:
    for line in handle:
        if not line.strip():
            continue

        barcode, sample = line.rstrip().split("\t")
        barcode_cb = barcode + "-1"
        cell = f"{barcode}.{sample}"
        folder = os.path.join(args.outdir, cell)
        os.makedirs(folder, exist_ok=True)

        with open(os.path.join(folder, "metadata.yaml"), "w") as meta:
            meta.write(f"barcode: {barcode}\n")
            meta.write(f"sample: {sample}\n")

        writers[barcode_cb] = gzip.open(
            os.path.join(folder, "reads.fq.gz"),
            "wt",
        )
        written[barcode_cb] = 0

print(
    f"[extract] bam={args.bam} cells={len(writers)} threads={args.threads}",
    flush=True,
)
t0 = time.time()

bam = pysam.AlignmentFile(args.bam, "rb", threads=args.threads)

n_seen = 0
n_kept = 0
for read in bam.fetch(until_eof=True):
    n_seen += 1
    if args.progress_every and n_seen % args.progress_every == 0:
        elapsed = time.time() - t0
        rate = n_seen / elapsed if elapsed > 0 else 0
        print(
            f"[extract] scanned={n_seen:,} kept={n_kept:,} "
            f"elapsed={elapsed:.0f}s rate={rate:,.0f} reads/s",
            flush=True,
        )

    if read.is_unmapped:
        continue

    try:
        cb = read.get_tag("CB")
    except KeyError:
        continue

    writer = writers.get(cb)
    if writer is None:
        continue

    writer.write(
        f"@{read.query_name}\n"
        f"{read.query_sequence}\n"
        f"+\n"
        f"{read.qual}\n"
    )
    written[cb] += 1
    n_kept += 1

bam.close()

for handle in writers.values():
    handle.close()

elapsed = time.time() - t0
print(
    f"[extract] done scanned={n_seen:,} kept={n_kept:,} "
    f"elapsed={elapsed:.0f}s cells="
    + ",".join(f"{cb}:{n}" for cb, n in written.items()),
    flush=True,
)
sys.exit(0)
