#!/usr/bin/env bash
# Build a pseudo-GTF from a transcriptome FASTA (wraps bin/create_gtf.sh).
# Usage: run_create_gtf.sh INPUT.fasta OUTPUT.gtf
set -euo pipefail

infile="${1:?}"
outfile="${2:?}"
root="$(cd "$(dirname "$0")/../.." && pwd)"

if ! command -v seqkit >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1; then
        module load seqkit || true
    fi
fi
if ! command -v seqkit >/dev/null 2>&1; then
    echo "ERROR: seqkit not found" >&2
    exit 1
fi

workdir="$(mktemp -d "${TMPDIR:-/tmp}/create_gtf.XXXXXX")"
trap 'rm -rf "$workdir"' EXIT

# create_gtf.sh writes 01_/02_/… relative to CWD using the literal FASTA arg name.
base="genes.final.fasta"
cp -f "$infile" "$workdir/$base"
(
    cd "$workdir"
    bash "$root/bin/create_gtf.sh" "$base"
    # Legacy strand fix from run_pipeline.sh
    sed -i 's/\.\t1\t\./.\t+\t./g' "04_${base}.gene_intervals.gtf" 2>/dev/null \
        || sed -i.bak 's/\.\t1\t\./.\t+\t./g' "04_${base}.gene_intervals.gtf"
)

mkdir -p "$(dirname "$outfile")"
cp -f "$workdir/04_${base}.gene_intervals.gtf" "$outfile"
echo "[gtf] wrote $outfile"
