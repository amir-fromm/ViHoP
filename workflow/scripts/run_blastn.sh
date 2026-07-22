#!/usr/bin/env bash
# BLAST 18S sequences against a nucleotide DB.
# Usage: run_blastn.sh QUERY.fasta DB OUT.tsv EVALUE THREADS
set -euo pipefail

query="${1:?}"
db="${2:?}"
out="${3:?}"
evalue="${4:-1e-10}"
threads="${5:-8}"

if ! command -v blastn >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1; then
        module purge || true
        module load BLAST/2.11.0-Linux_x86_64 || module load BLAST || true
    fi
fi

if ! command -v blastn >/dev/null 2>&1; then
    echo "ERROR: blastn not found" >&2
    exit 1
fi

mkdir -p "$(dirname "$out")"

if [[ ! -s "$query" ]]; then
    echo "[blastn] empty query; writing empty $out"
    : > "$out"
    exit 0
fi

echo "[blastn] query=$query db=$db threads=$threads"
blastn -outfmt 6 -evalue "$evalue" -query "$query" -db "$db" \
    -num_threads "$threads" -out "$out"
echo "[blastn] wrote $out ($(wc -l < "$out") hits)"
