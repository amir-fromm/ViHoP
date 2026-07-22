#!/usr/bin/env bash
# Per-cell 18S LAST search on extracted reads (legacy lastal_array_18s_run.sh).
# Usage: run_lastal_18s.sh READS.fq.gz OUT.last DB IDENTITY LENGTH THREADS MEM_MB
set -euo pipefail

reads="${1:?}"
out="${2:?}"
db="${3:?}"
identity="${4:-99}"
length="${5:-54}"
threads="${6:-2}"
mem_mb="${7:-4000}"

if ! command -v lastal >/dev/null 2>&1 || ! command -v seqtk >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1; then
        module purge || true
        module load seqtk || true
        # last often comes from conda; try module if present
        module load last || true
    fi
fi

if ! command -v lastal >/dev/null 2>&1; then
    echo "ERROR: lastal not found" >&2
    exit 1
fi
if ! command -v seqtk >/dev/null 2>&1; then
    echo "ERROR: seqtk not found" >&2
    exit 1
fi

mkdir -p "$(dirname "$out")"

# lastal -i expects bytes; convert MB → rough GB floor of at least 1.
mem_gb=$(( mem_mb / 1000 ))
if [[ "$mem_gb" -lt 1 ]]; then
    mem_gb=1
fi

echo "[lastal] reads=$reads db=$db identity=$identity length=$length"
seqtk seq -a "$reads" | \
    lastal -m 100 -v -f BlastTab -P0 -i"${mem_gb}G" -u 2 -Q 0 "$db" | \
    awk -v id="$identity" -v len="$length" '$0 !~ /^#/ && $3 >= id && $4 >= len' \
    > "$out"

echo "[lastal] wrote $out ($(wc -l < "$out") hits)"
