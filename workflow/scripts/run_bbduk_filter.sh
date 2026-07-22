#!/usr/bin/env bash
# Filter low-complexity transcripts with BBDuk (legacy run_pipeline.sh step).
# Usage: run_bbduk_filter.sh INPUT.fasta OUTPUT.fasta OUTM.fq ENTROPY
set -euo pipefail

infile="${1:?input fasta required}"
outfile="${2:?output fasta required}"
outm="${3:?low-complexity output required}"
entropy="${4:-0.7}"

if ! command -v bbduk.sh >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1; then
        module purge || true
        module load BBMap || true
    fi
fi

if ! command -v bbduk.sh >/dev/null 2>&1; then
    echo "ERROR: bbduk.sh not found in PATH or modules" >&2
    exit 1
fi

mkdir -p "$(dirname "$outfile")" "$(dirname "$outm")"

echo "[bbduk] in=$infile entropy=$entropy"
bbduk.sh in="$infile" out="$outfile" outm="$outm" entropy="$entropy"
echo "[bbduk] wrote $outfile"
