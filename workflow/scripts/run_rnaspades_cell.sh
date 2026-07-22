#!/usr/bin/env bash
# Option A: assemble one cell from a single reads.fq.gz via RNA-SPAdes --s1.
# Usage: run_rnaspades_cell.sh READS_FQ_GZ OUTPUT_TRANSCRIPTS THREADS
set -euo pipefail

reads="${1:?reads.fq.gz required}"
output="${2:?output transcripts.fasta required}"
threads="${3:-1}"

if [[ ! -s "$reads" ]]; then
    echo "ERROR: missing or empty reads file: $reads" >&2
    exit 1
fi

# Prefer PATH (conda); fall back to Wexac module as in the legacy pipeline.
if ! command -v rnaspades.py >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1; then
        module purge || true
        module load SPAdes/3.15.3-GCC-11.2.0 || module load SPAdes || true
    fi
fi

if ! command -v rnaspades.py >/dev/null 2>&1; then
    echo "ERROR: rnaspades.py not found in PATH or modules" >&2
    exit 1
fi

outdir="$(dirname "$output")"
mkdir -p "$outdir"

workdir="$(mktemp -d "${TMPDIR:-/tmp}/vihop_spades.XXXXXX")"
cleanup() { rm -rf "$workdir"; }
trap cleanup EXIT

echo "[assemble] reads=$reads threads=$threads workdir=$workdir"
rnaspades.py -k auto -t "$threads" -o "$workdir" --s1 "$reads"

if [[ ! -s "$workdir/transcripts.fasta" ]]; then
    echo "ERROR: SPAdes did not produce transcripts.fasta" >&2
    ls -la "$workdir" >&2 || true
    exit 1
fi

cp "$workdir/transcripts.fasta" "$output"
echo "[assemble] wrote $output ($(wc -c < "$output") bytes)"
