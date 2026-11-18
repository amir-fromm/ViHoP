#!/usr/bin/env bash
# MMseqs local-alignment clustering (wraps logic of bin/mmseqs_cluster.sh).
# Usage: run_mmseqs_dedup.sh INPUT.fasta OUTPUT.fasta IDENTITY THREADS \
#          [COV] [COV_MODE] [MIN_ALN_LEN] [WORKDIR]
set -euo pipefail

infile="${1:?}"
outfile="${2:?}"
identity="${3:-0.95}"
threads="${4:-32}"
cov="${5:-0.30}"
cov_mode="${6:-1}"
min_aln_len="${7:-200}"
workdir="${8:-}"

if ! command -v mmseqs >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1; then
        module purge || true
        module load MMseqs2 || module load mmseqs2 || module load mmseqs || true
    fi
fi
if ! command -v mmseqs >/dev/null 2>&1; then
    echo "ERROR: mmseqs not found" >&2
    exit 1
fi

mkdir -p "$(dirname "$outfile")"
if [[ -z "$workdir" ]]; then
    workdir="$(dirname "$outfile")/mmseqs_work"
fi
tmpdir="${TMPDIR:-/tmp}/mmseqs_$$"
mkdir -p "$workdir" "$tmpdir"
trap 'rm -rf "$tmpdir"' EXIT

echo "[mmseqs] in=$infile identity=$identity cov=$cov cov_mode=$cov_mode min_aln_len=$min_aln_len threads=$threads"

mmseqs createdb "$infile" "$workdir/db"
mmseqs cluster \
    "$workdir/db" "$workdir/clu" "$tmpdir" \
    --min-seq-id "$identity" \
    -c "$cov" --cov-mode "$cov_mode" \
    --min-aln-len "$min_aln_len" \
    --threads "$threads"
mmseqs createtsv "$workdir/db" "$workdir/db" "$workdir/clu" "$workdir/clusters.tsv"
mmseqs result2repseq "$workdir/db" "$workdir/clu" "$workdir/repdb" --threads "$threads"
mmseqs result2flat "$workdir/db" "$workdir/db" "$workdir/repdb" "$workdir/dedup_reps.fa" --use-fasta-header

cp -f "$workdir/dedup_reps.fa" "$outfile"
echo "[mmseqs] wrote $outfile ($(grep -c '^>' "$outfile" || true) seqs)"
echo "[mmseqs] cluster table: $workdir/clusters.tsv"
