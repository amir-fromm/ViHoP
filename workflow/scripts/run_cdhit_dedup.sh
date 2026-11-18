#!/usr/bin/env bash
# CD-HIT transcriptome deduplication.
# Usage: run_cdhit_dedup.sh INPUT.fasta OUTPUT.fasta IDENTITY THREADS
set -euo pipefail

infile="${1:?}"
outfile="${2:?}"
identity="${3:-0.95}"
threads="${4:-16}"

if ! command -v cd-hit-est >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1; then
        module purge || true
        module load CD-HIT/4.6.6-GCC-10.2.0 || module load CD-HIT || true
    fi
fi
if ! command -v cd-hit-est >/dev/null 2>&1; then
    echo "ERROR: cd-hit-est not found" >&2
    exit 1
fi

mkdir -p "$(dirname "$outfile")"
echo "[cdhit] in=$infile identity=$identity threads=$threads"
cd-hit-est -c "$identity" -T "$threads" -M 0 -i "$infile" -o "$outfile"
echo "[cdhit] wrote $outfile ($(grep -c '^>' "$outfile" || true) seqs)"
