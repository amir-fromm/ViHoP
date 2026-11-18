#!/usr/bin/env bash
# Filter low-complexity transcripts with BBDuk (legacy run_pipeline.sh step).
# Usage: run_bbduk_filter.sh INPUT.fasta OUTPUT.fasta OUTM.fq ENTROPY [MEM_MB]
set -euo pipefail

infile="${1:?input fasta required}"
outfile="${2:?output fasta required}"
outm="${3:?low-complexity output required}"
entropy="${4:-0.7}"
mem_mb="${5:-4000}"

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

# BBMap otherwise sizes the heap from the whole node. On a 325 GB node that
# request exceeded the LSF reservation and the job was killed (TERM_MEMLIMIT).
if ! [[ "$mem_mb" =~ ^[0-9]+$ ]]; then
    echo "ERROR: memory MB must be an integer, got: $mem_mb" >&2
    exit 1
fi
heap_mb=$(( mem_mb - 1024 ))
if (( heap_mb < 512 )); then
    heap_mb=512
fi

echo "[bbduk] in=$infile entropy=$entropy heap=${heap_mb}m"
bbduk.sh -Xmx"${heap_mb}m" in="$infile" out="$outfile" outm="$outm" entropy="$entropy"
echo "[bbduk] wrote $outfile"
